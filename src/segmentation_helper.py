#!/usr/bin/env python3

import warnings
warnings.filterwarnings('ignore')
import tensorflow as tf
from tensorflow.keras import backend as K
from tensorflow.keras.preprocessing.sequence import pad_sequences
from tensorflow.keras.models import load_model
from tensorflow.keras.metrics import MeanIoU


# Padding
def compute_padding_length(reads_lengths, nb_pooling=5, pooling_factor=2):
    # bugfix when reads_lengths is a vector of length 1 in R
    if not isinstance(reads_lengths, list):
        reads_lengths = [reads_lengths]
    # compute padding length
    div_length = (pooling_factor**nb_pooling)
    max_read_length = max(reads_lengths)
    max_pad_length = (max_read_length//div_length+1)*div_length
    return max_pad_length

def pad_reads(reads_list, max_pad_length):
    # transform list to 3D array & add padding
    reads_array = pad_sequences(reads_list, dtype="float", padding="post", maxlen=max_pad_length)
    return reads_array

# IOU metric - https://gist.github.com/ilmonteux/8340df952722f3a1030a7d937e701b5a
# DUPLICATED
def seg_metrics(y_true, y_pred, metric_name, metric_type='standard', drop_last = True, mean_per_class=False, verbose=False):
    """ 
    Compute mean metrics of two segmentation masks, via Keras.
    
    IoU(A,B) = |A & B| / (| A U B|)
    Dice(A,B) = 2*|A & B| / (|A| + |B|)
    
    Args:
        y_true: true masks, one-hot encoded.
        y_pred: predicted masks, either softmax outputs, or one-hot encoded.
        metric_name: metric to be computed, either 'iou' or 'dice'.
        metric_type: one of 'standard' (default), 'soft', 'naive'.
          In the standard version, y_pred is one-hot encoded and the mean
          is taken only over classes that are present (in y_true or y_pred).
          The 'soft' version of the metrics are computed without one-hot 
          encoding y_pred.
          The 'naive' version return mean metrics where absent classes contribute
          to the class mean as 1.0 (instead of being dropped from the mean).
        drop_last = True: boolean flag to drop last class (usually reserved
          for background class in semantic segmentation)
        mean_per_class = False: return mean along batch axis for each class.
        verbose = False: print intermediate results such as intersection, union
          (as number of pixels).
    Returns:
        IoU/Dice of y_true and y_pred, as a float, unless mean_per_class == True
          in which case it returns the per-class metric, averaged over the batch.
    
    Inputs are B*W*H*N tensors, with
        B = batch size,
        W = width,
        H = height,
        N = number of classes
    """
    
    flag_soft = (metric_type == 'soft')
    flag_naive_mean = (metric_type == 'naive')
    
    # always assume one or more classes
    num_classes = K.shape(y_true)[-1]
        
    if not flag_soft:
        # get one-hot encoded masks from y_pred (true masks should already be one-hot)
        y_pred = K.one_hot(K.argmax(y_pred), num_classes)
        y_true = K.one_hot(K.argmax(y_true), num_classes)

    # if already one-hot, could have skipped above command
    # keras uses float32 instead of float64, would give error down (but numpy arrays or keras.to_categorical gives float64)
    y_true = K.cast(y_true, 'float32')
    y_pred = K.cast(y_pred, 'float32')

    # intersection and union shapes are batch_size * n_classes (values = area in pixels)
    axes = (1) # W,H axes of each image
    intersection = K.sum(K.abs(y_true * y_pred), axis=axes)
    mask_sum = K.sum(K.abs(y_true), axis=axes) + K.sum(K.abs(y_pred), axis=axes)
    union = mask_sum - intersection # or, np.logical_or(y_pred, y_true) for one-hot

    smooth = .001
    iou = (intersection + smooth) / (union + smooth)
    dice = 2 * (intersection + smooth)/(mask_sum + smooth)

    metric = {'iou': iou, 'dice': dice}[metric_name]

    # define mask to be 0 when no pixels are present in either y_true or y_pred, 1 otherwise
    mask =  K.cast(K.not_equal(union, 0), 'float32')
    
    if drop_last:
        metric = metric[:,:-1]
        mask = mask[:,:-1]
    
    if verbose:
        print('intersection, union')
        print(K.eval(intersection), K.eval(union))
        print(K.eval(intersection/union))
    
    # return mean metrics: remaining axes are (batch, classes)
    if flag_naive_mean:
        return K.mean(metric)

    # take mean only over non-absent classes
    class_count = K.sum(mask, axis=0)
    non_zero = tf.greater(class_count, 0)
    non_zero_sum = tf.boolean_mask(K.sum(metric * mask, axis=0), non_zero)
    non_zero_count = tf.boolean_mask(class_count, non_zero)
    
    if verbose:
        print('Counts of inputs with class present, metrics for non-absent classes')
        print(K.eval(class_count), K.eval(non_zero_sum / non_zero_count))
        
    return K.mean(non_zero_sum / non_zero_count)

def mean_iou(y_true, y_pred, **kwargs):
    """
    Compute mean Intersection over Union of two segmentation masks, via Keras.

    Calls metrics_k(y_true, y_pred, metric_name='iou'), see there for allowed kwargs.
    """
    return seg_metrics(y_true, y_pred, metric_name='iou', metric_type="standard", drop_last=False, **kwargs)

# Define IOU loss
# https://github.com/karolzak/keras-unet/blob/master/keras_unet/losses.py
def jaccard_distance(y_true, y_pred, smooth=100):
    """Jaccard distance for semantic segmentation.
    Also known as the intersection-over-union loss.
    This loss is useful when you have unbalanced numbers of pixels within an image
    because it gives all classes equal weight. However, it is not the defacto
    standard for image segmentation.
    For example, assume you are trying to predict if
    each pixel is cat, dog, or background.
    You have 80% background pixels, 10% dog, and 10% cat.
    If the model predicts 100% background
    should it be be 80% right (as with categorical cross entropy)
    or 30% (with this loss)?
    The loss has been modified to have a smooth gradient as it converges on zero.
    This has been shifted so it converges on 0 and is smoothed to avoid exploding
    or disappearing gradient.
    Jaccard = (|X & Y|)/ (|X|+ |Y| - |X & Y|)
            = sum(|A*B|)/(sum(|A|)+sum(|B|)-sum(|A*B|))
    # Arguments
        y_true: The ground truth tensor.
        y_pred: The predicted tensor
        smooth: Smoothing factor. Default is 100.
    # Returns
        The Jaccard distance between the two tensors.
    # References
        - [What is a good evaluation measure for semantic segmentation?](
           http://www.bmva.org/bmvc/2013/Papers/paper0032/paper0032.pdf)
    """
    intersection = K.sum(K.abs(y_true * y_pred), axis=-1)
    sum_ = K.sum(K.abs(y_true) + K.abs(y_pred), axis=-1)
    jac = (intersection + smooth) / (sum_ - intersection + smooth)
    return (1 - jac) * smooth

def rle_segments(mask, class_value):
    """
    Find segments of consecutive indices where the mask equals the class_value using RLE.
    
    Args:
        mask: Array of class predictions.
        class_value: The class value (0 or 1) for which to find segments.
    
    Returns:
        A list of tuples (start_idx, end_idx) for each segment.
    """
    # Create a boolean mask where the mask is equal to the class_value
    is_class = (mask == class_value)
    
    # Calculate the differences between adjacent values in the boolean mask
    diff = np.diff(is_class.astype(int))
    
    # The start of a segment occurs where the diff is 1, and the end occurs where diff is -1
    start_indices = np.where(diff == 1)[0] + 1  # shift by 1 due to diff
    end_indices = np.where(diff == -1)[0]
    
    # Handle the case where the mask starts with the class_value
    if is_class[0]:
        start_indices = np.insert(start_indices, 0, 0)
    
    # Handle the case where the mask ends with the class_value
    if is_class[-1]:
        end_indices = np.append(end_indices, len(mask) - 1)
    
    return list(zip(start_indices, end_indices))

def compute_x0_distance(y_true, y_pred):
    y_true_classes = np.argmax(y_true, axis=1)
    y_pred_classes = np.argmax(y_pred, axis=1)
    
    distances = []
    positions = []
    for class_id in [0, 1]:  # Only process class 0 and 1
        true_segments = rle_segments(y_true_classes, class_id)
        pred_segments = rle_segments(y_pred_classes, class_id)
        
        # Compute distances between each segment in y_true and the corresponding segment in y_pred
        for pred_segment in pred_segments:
            for true_segment in true_segments:
                if pred_segment[0] <= true_segment[1] and pred_segment[1] >= true_segment[0]:
                    x0_true = true_segment[0] if class_id == 0 else true_segment[1]
                    x0_pred = pred_segment[0] if class_id == 0 else pred_segment[1]
                    distances.append(abs(x0_true - x0_pred))
                    positions.append((x0_true + x0_pred)/2)
    if not distances:
        # return 0.0, 0.0  # If no valid distances are found
        return [], []  # If no valid distances are found
    
    return distances, positions

def compute_x0_distance_batch(y_true_batch, y_pred_batch):
    all_distances = []
    all_positions = []
    
    for i in range(len(y_true_batch)):  # Loop through each sample in the batch
        distances, positions = compute_x0_distance(y_true_batch[i], y_pred_batch[i])
        all_distances.append(distances)
        all_positions.append(positions)
        # all_distances.extend(distances)
        # all_positions.extend(positions)

    if not all_distances:
        # return 0.0, 0.0 # If no valid distances are found
        return [] # If no valid distances are found
    
    return all_distances

class MeanX0Distance(tf.keras.metrics.Metric):
    def __init__(self, name="mean_x0_distance", **kwargs):
        super(MeanX0Distance, self).__init__(name=name, **kwargs)
        self.total_distance = self.add_weight(name="total", initializer="zeros")
        self.count = self.add_weight(name="count", initializer="zeros")
    
    def update_state(self, y_true, y_pred, sample_weight=None):
        distances = tf.py_function(
            func=compute_x0_distance_batch, inp=[y_true, y_pred], Tout=tf.float32
        )
        # if distances:
        self.total_distance.assign_add(tf.reduce_sum(distances))
        self.count.assign_add(tf.cast(tf.size(distances), tf.float32))
    
    def result(self):
        return self.total_distance / self.count
    
    def reset_state(self):
        self.total_distance.assign(0.0)
        self.count.assign(0.0)
# DUPLICATED

# Prediction
def predict(reads_array, path_to_model, dependencies={"mean_iou": mean_iou, 'MeanX0Distance': MeanX0Distance, 'MeanIoU': MeanIoU}):
    neuralnet = load_model(path_to_model, custom_objects=dependencies)
    predictions = neuralnet.predict(reads_array)
    
    return predictions


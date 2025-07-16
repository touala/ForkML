# U-Net
# Adapted from this script : 
# https://github.com/karolzak/keras-unet/blob/master/keras_unet/models/custom_unet.py

from tensorflow.keras.models import Model
from tensorflow.keras.layers import (BatchNormalization, Conv1D, MaxPooling1D, Dropout, SpatialDropout1D, UpSampling1D, Input, concatenate, multiply, add, Activation, Masking)

def conv1d_block(
    inputs,
    use_batch_norm=True,
    dropout=0.3,
    dropout_type="spatial",
    filters=16,
    kernel_size=10,
    activation="relu",
    kernel_initializer="he_normal",
    padding="same"
    ):
    
    if dropout_type == "spatial":
        DO = SpatialDropout1D
    elif dropout_type == "standard":
        DO = Dropout
    else:
        raise ValueError(
            f"dropout_type must be one of ['spatial', 'standard'], got {dropout_type}"
        )
    
    c = Conv1D(
        filters,
        kernel_size,
        activation=activation,
        kernel_initializer=kernel_initializer,
        padding=padding,
        use_bias=not use_batch_norm,
    )(inputs)
    if use_batch_norm:
        c = BatchNormalization()(c)
    if dropout > 0.0:
        c = DO(dropout)(c)
    c = Conv1D(
        filters,
        kernel_size,
        activation=activation,
        kernel_initializer=kernel_initializer,
        padding=padding,
        use_bias=not use_batch_norm,
    )(c)
    if use_batch_norm:
        c = BatchNormalization()(c)
    return c

def custom_unet(
    input_shape, num_classes=2, activation="relu", use_batch_norm=True,
    dropout=0.3, dropout_change_per_layer=0.0, dropout_type="spatial",
    use_dropout_on_upsampling=False, filters=16, num_layers=4,
    output_activation="sigmoid", kernel_size=10, include_masking="default"
    ):
    
    """
    Customisable UNet architecture (Ronneberger et al. 2015 [1]).
    Arguments:
    input_shape: 3D Tensor of shape (x, y, num_channels)
    num_classes (int): Unique classes in the output mask. Should be set to 1 for binary segmentation
    activation (str): A keras.activations.Activation to use. ReLu by default.
    use_batch_norm (bool): Whether to use Batch Normalisation across the channel axis between convolutional layers
    dropout (float between 0. and 1.): Amount of dropout after the initial convolutional block. Set to 0. to turn Dropout off
    dropout_change_per_layer (float between 0. and 1.): Factor to add to the Dropout after each convolutional block
    dropout_type (one of "spatial" or "standard"): Type of Dropout to apply. Spatial is recommended for CNNs [2]
    use_dropout_on_upsampling (bool): Whether to use dropout in the decoder part of the network
    filters (int): Convolutional filters in the initial convolutional block. Will be doubled every block
    num_layers (int): Number of total layers in the encoder not including the bottleneck layer
    output_activation (str): A keras.activations.Activation to use. Sigmoid by default for binary segmentation
    Returns:
    model (keras.models.Model): The built U-Net
    Raises:
    ValueError: If dropout_type is not one of "spatial" or "standard"
    [1]: https://arxiv.org/abs/1505.04597
    [2]: https://arxiv.org/pdf/1411.4280.pdf
    [3]: https://arxiv.org/abs/1804.03999
    """
    
    # Build U-Net model
    inputs = Input(input_shape)
    
    if include_masking == "clean":
        x = Masking(mask_value=-1)(inputs)  # Masking padding values
    elif include_masking == "rough":
        x = Masking(mask_value=0.0)(inputs)  # Masking padding values    
    else:
        x = inputs
    
    # Down layers
    down_layers = []
    for l in range(num_layers):
        # print(f'{l} {filters}')
        x = conv1d_block(
            inputs=x,
            filters=filters,
            use_batch_norm=use_batch_norm,
            dropout=dropout,
            dropout_type=dropout_type,
            activation=activation,
            kernel_size=kernel_size
        )
        down_layers.append(x)
        x = MaxPooling1D(2)(x)
        dropout += dropout_change_per_layer
        filters = filters * 2  # double the number of filters with each layer
    
    # print(f'{l} {filters}')
    x = conv1d_block(
        inputs=x,
        filters=filters,
        use_batch_norm=use_batch_norm,
        dropout=dropout,
        dropout_type=dropout_type,
        activation=activation,
        kernel_size=kernel_size
    )
    
    # Up layers
    if not use_dropout_on_upsampling:
        dropout = 0.0
        dropout_change_per_layer = 0.0
    
    for conv in reversed(down_layers):
        filters //= 2  # decreasing number of filters with each layer
        dropout -= dropout_change_per_layer
        x = UpSampling1D(2)(x)
        x = concatenate([x, conv])
        x = conv1d_block(
            inputs=x,
            filters=filters,
            use_batch_norm=use_batch_norm,
            dropout=dropout,
            dropout_type=dropout_type,
            activation=activation,
            kernel_size=kernel_size
        )
        # print(f'{conv} {filters}')
    
    outputs = Conv1D(num_classes, 1, activation=output_activation)(x)
    
    # Return model
    model = Model(inputs=[inputs], outputs=[outputs])
    return model

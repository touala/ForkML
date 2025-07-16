#!/bin/bash

ln -s $CONDA_PREFIX/lib/python3.10/site-packages/tensorrt/libnvinfer.so.8 $CONDA_PREFIX/lib/python3.10/site-packages/tensorrt/libnvinfer.so.7
ln -s $CONDA_PREFIX/lib/python3.10/site-packages/tensorrt/libnvinfer_plugin.so.8 $CONDA_PREFIX/lib/python3.10/site-packages/tensorrt/libnvinfer_plugin.so.7

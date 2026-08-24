#!/bin/bash

srun \
  --partition=gpu \
  --gres=gpu:a100_40gb:1 \
  --cpus-per-task=4 \
  --mem=32G \
  --time=02:00:00 \
  --pty bash

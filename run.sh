#!/bin/bash

nohup Rscript pipeline_patch.R $1  > salida.out 2>&1 &

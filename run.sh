#!/bin/bash

nohup Rscript pipeline_mod.R $1  > salida.out 2>&1 &

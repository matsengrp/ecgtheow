#!/bin/sh

set -eux

for i in data inferences runs; do
    rsync -avz stoat:/home/matsengrp/working/matsen/ecgtheow/$i/ $i
done


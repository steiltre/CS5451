#!/bin/bash

if [ "$#" -ne 1 ]; then
  echo "usage: ./plot.sh <data file>";
  echo "  NOTE: 'clusters.txt' must be in working directory.";
  exit 1;
fi

octave-cli --eval "plot_clusters2D('$1', 'clusters.txt', 'centroids.txt')"


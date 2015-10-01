#! /bin/bash
# The PYTHONPATH must point to the ilastik source folder and the installed pgmlink library

SCRIPT_DIR=$(cd `dirname $0` && pwd)

python ${SCRIPT_DIR}/../workflow/workflow.py ${SCRIPT_DIR}/testData/Sequence000-003_probabilities.h5 -r ${SCRIPT_DIR}/testData/Sequence000-003.h5 --prediction-channel=1 --normalize-predictions --sigma=0.3 --background-threshold=128 --t-axis=0 --ch-axis=3 --percentile=0.8 --progressive=0.9 --num-levels=4 --detectionClassifier=${SCRIPT_DIR}/testData/ilastik-object-detection.ilp --seedRadius=1 --classifierFile=${SCRIPT_DIR}/testData/classifiers.h5 --output=${SCRIPT_DIR}/result.h5 --outputOversegmentation=${SCRIPT_DIR}/testData/oversegmentation.h5 --oversegmentationOnly

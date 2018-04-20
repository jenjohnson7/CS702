#!/bin/bash

readonly JAR_URL='https://linqs-data.soe.ucsc.edu/maven/repositories/psl-releases/org/linqs/psl-cli/CANARY-2.1.1/psl-cli-CANARY-2.1.1.jar'
readonly JAR_PATH='psl-cli-CANARY-2.1.1.jar'

function main() {
   trap exit SIGINT

   # Make sure we can run PSL.
   check_requirements
   fetch_file "${JAR_URL}" "${JAR_PATH}" 'psl-jar'

   # Run PSL

  #  AN
   runWeightLearning "AN_sum" "AN"
  #  runEvaluation "AN_sum" "AN" "AN_sum_results"
  #
  #  runWeightLearning "AN_no_sum" "AN"
  #  runEvaluation "AN_no_sum" "AN" "AN_no_sum_results"
  # 
  # #  No AN
  #
  #  runWeightLearning "no_AN_sum" "no_AN"
  #  runEvaluation "no_AN_sum" "no_AN" "no_AN_sum_results"
  #
  #  runWeightLearning "no_AN_no_sum" "no_AN"
  #  runEvaluation "no_AN_no_sum" "no_AN" "no_AN_no_sum_results"

 }

function runWeightLearning() {
   echo "Running PSL Weight Learning"

   local model=$1
   local data=$2

   java -jar "${JAR_PATH}" -learn -model "${model}.psl" -data "${data}.data"
   if [[ "$?" -ne 0 ]]; then
      echo 'ERROR: Failed to run weight learning'
      exit 60
   fi
}

function runEvaluation() {
   echo "Running PSL Inference"

   local model=$1
   local data=$2
   local results=$3

   java -jar "${JAR_PATH}" -infer -model "${model}-learned.psl" -data "${data}.data" -output "../output_data/"${results}
   if [[ "$?" -ne 0 ]]; then
      echo 'ERROR: Failed to run infernce'
      exit 70
   fi
}

function check_requirements() {
   local hasWget
   local hasCurl

   type wget > /dev/null 2> /dev/null
   hasWget=$?

   type curl > /dev/null 2> /dev/null
   hasCurl=$?

   if [[ "${hasWget}" -ne 0 ]] && [[ "${hasCurl}" -ne 0 ]]; then
      echo 'ERROR: wget or curl required to download dataset'
      exit 10
   fi

   type java > /dev/null 2> /dev/null
   if [[ "$?" -ne 0 ]]; then
      echo 'ERROR: java required to run project'
      exit 13
   fi
}

function get_fetch_command() {
   type curl > /dev/null 2> /dev/null
   if [[ "$?" -eq 0 ]]; then
      echo "curl -o"
      return
   fi

   type wget > /dev/null 2> /dev/null
   if [[ "$?" -eq 0 ]]; then
      echo "wget -O"
      return
   fi

   echo 'ERROR: wget or curl not found'
   exit 20
}

function fetch_file() {
   local url=$1
   local path=$2
   local name=$3

   if [[ -e "${path}" ]]; then
      echo "${name} file found cached, skipping download."
      return
   fi

   echo "Downloading ${name} file with command: $FETCH_COMMAND"
   `get_fetch_command` "${path}" "${url}"
   if [[ "$?" -ne 0 ]]; then
      echo "ERROR: Failed to download ${name} file"
      exit 30
   fi
}

main "$@"

#!/bin/bash

available_targets=(advx advy difx dify pres)

target=${1}

# check if target is an expected value
expected=false
for t in ${available_targets[@]}; do
  if [[ ${t} == ${target} ]]; then
    expected=true
    break
  fi
done
if [ ${expected} = false ]; then
  echo unexpected argument: ${target}
  exit 1
fi

resols=(8 16 32 64 128)

: > ${target}.dat
for resol in ${resols[@]}; do
  cc \
    -DTEST \
    -std=c99 -Wall -Wextra \
    -I../../../../include \
    ../../../memory.c \
    ../../../array.c \
    ../../../domain.c \
    ../test_util.c \
    ./test_util.c \
    ${target}.c \
    -o a.out \
    -lm \
    && \
    ./a.out ${resol} >> ${target}.dat
done


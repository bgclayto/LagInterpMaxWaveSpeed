#!/bin/bash

./lambda > tests/output.txt
diff tests/output.txt tests/test2.txt

if [ $? -eq 0 ]
then
    exit 0
else
    exit 1
fi


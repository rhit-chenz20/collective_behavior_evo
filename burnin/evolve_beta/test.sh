#!/bin/bash
c=0
for rep in {1..20}; do  
    ((c++))
    echo "count is $c"
    if (( c > 10)); then
        sleep 2
        c=0
    fi
done
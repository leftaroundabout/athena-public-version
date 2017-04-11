#!/bin/bash
./$(grep config run/athinput | sed 's/ =/.py/')
make clean
make

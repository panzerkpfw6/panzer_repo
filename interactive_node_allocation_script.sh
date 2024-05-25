#!/bin/bash
salloc -n 1 -A k1205 --partition=workq --time=12:00:00 --hint=nomultithread --exclusive
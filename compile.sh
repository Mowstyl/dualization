#!/bin/bash

cd builddir/ && meson compile | tee ../compile.log ; cd ..

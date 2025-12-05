#!/bin/bash

if [[ ! -d ./builddir ]]; then
    meson builddir
fi

cd builddir/ && meson compile | tee ../compile.log ; cd ..

{ pkgs ? import <nixpkgs> {} }:

pkgs.mkShell {
  buildInputs = [
    pkgs.python39
    pkgs.python39Packages.numpy
    pkgs.python39Packages.matplotlib
    pkgs.python39Packages.numba
    pkgs.python39Packages.tkinter
    pkgs.poetry
  ];

  shellHook = ''
    export VIRTUAL_ENV=$PWD/.venv
    export PATH=$VIRTUAL_ENV/bin:$PATH
    if [ ! -d "$VIRTUAL_ENV" ]; then
      python -m venv $VIRTUAL_ENV
    fi
    source $VIRTUAL_ENV/bin/activate
  '';
}


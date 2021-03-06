#!/bin/bash

GIT_LFS_VERSION=v0.6.0

echo "cloning"
git clone https://github.com/github/git-lfs.git

echo "cd into git-lfs"
cd git-lfs

echo "checkout commit e8a5cfb85e27c1284a73a2eb8269bd5f5c4d0955"
git checkout $GIT_LFS_VERSION

echo "scripts/bootstrap"
script/bootstrap

echo "ls bin"
ls bin

echo "cd back down"
cd ..

echo "resetting travis remote"
git remote set-url origin "git@github.com:broadinstitute/hellbender.git"

echo "init"
git-lfs/bin/git-lfs init

echo "pull"
git-lfs/bin/git-lfs pull

echo "ls-files"
git-lfs/bin/git-lfs ls-files

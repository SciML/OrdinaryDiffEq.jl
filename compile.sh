#!/bin/bash

# Set the path to your Julia project
PROJECT_DIR="C:\Users\srira\Downloads\ode\OrdinaryDiffEq.jl"

# Change to the project directory
cd "$PROJECT_DIR" || { echo "Directory not found: $PROJECT_DIR"; exit 1; }

# Activate the project environment and precompile
julia --project=. -e 'using Pkg; Pkg.precompile()'

# Run your main Julia file (adjust the path as needed)
julia --project=. src/OrdinaryDiffEq.jl
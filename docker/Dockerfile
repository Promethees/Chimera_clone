FROM continuumio/miniconda3

# Install system dependencies for PyMOL (OpenGL and related libraries)
RUN apt-get update && apt-get install -y \
    libgl1-mesa-glx \
    libglu1-mesa \
    libglew-dev \
    libpng-dev \
    libfreetype6-dev \
    && rm -rf /var/lib/apt/lists/*

# Copy the environment.yml file
COPY environment.yml .

# Create the conda environment
RUN conda env create -f environment.yml

# Activate the conda environment by default
RUN echo "source activate chimera_clone" > ~/.bashrc

# Copy the project files
COPY . /app
WORKDIR /app

# Command to run the pipeline
CMD ["bash", "scripts/run_pipeline.sh"]
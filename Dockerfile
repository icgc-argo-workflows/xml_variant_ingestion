# Use a base image with Python pre-installed (assumed, since pip is used later)
FROM python:3.12.3-slim

# Set environment variable for consistent locale and encoding behavior
ENV LANG=C.UTF-8

# Set the working directory inside the container
# All subsequent commands will be run from this directory
WORKDIR /app

# Copy all files from the current directory (on host) to /app in the container
#COPY . /app

# Install the 'requests' package using pip without caching to reduce image size
#COPY requirements.txt /app/requirements.txt
#COPY requirements.txt requirements.txt
RUN apt-get update && apt-get install -y procps && rm -rf /var/lib/apt/lists/*
RUN pip install --no-cache-dir pandas biopython pyfaidx requests

# Set the default executable for the container
# This script will run first and usually wraps or prepares the environment
#ENTRYPOINT [ "/usr/local/env-execute" ]

# Provide default arguments to the ENTRYPOINT
# If no command is provided, this will effectively run: /usr/local/env-execute bash
CMD [ "bash" ]
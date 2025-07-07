# Start with Python 3.12 base image
FROM python:3.12-slim

# Set environment variables to prevent Python from buffering output
ENV PYTHONUNBUFFERED=1
ENV PYTHONDONTWRITEBYTECODE=1

RUN apt-get update && apt-get install -y \
    build-essential \
    python3-dev

RUN apt-get update && apt-get install -y libxrender1 libxext-dev git && rm -rf /var/lib/apt/lists/*
 
COPY * ./app/
WORKDIR /app/

RUN pip install --no-cache-dir -r requirements.txt 

# Expose the port the app runs on
EXPOSE 8000
# Command to run the application
CMD ["streamlit", "run", "main.py", "--server.address", "0.0.0.0", "--server.port", "8671"]

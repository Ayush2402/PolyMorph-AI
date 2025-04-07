#!/bin/bash

# Create Python virtual environment
python -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate

# Install Python dependencies
pip install -r requirements.txt

# Install frontend dependencies
cd frontend
npm install

# Start backend server (in background)
cd ../backend
uvicorn main:app --reload &

# Start frontend development server
cd ../frontend
npm start 
@echo off

REM Create Python virtual environment
python -m venv venv
call venv\Scripts\activate

REM Install Python dependencies
pip install -r requirements.txt

REM Install frontend dependencies
cd frontend
call npm install

REM Start backend server (in new window)
start cmd /k "cd backend && ..\venv\Scripts\uvicorn main:app --reload"

REM Start frontend development server
cd ../frontend
npm start 
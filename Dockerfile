FROM python:3.11-slim

WORKDIR /app

RUN apt-get update && apt-get install -y \
    build-essential \
    curl \
    software-properties-common \
    git \
    && rm -rf /var/lib/apt/lists/*

COPY requirements.txt requirements.txt
COPY app.py app.py
COPY pages/ pages/
COPY utils/ utils/
COPY resources/ resources/

RUN python3 -m pip install --upgrade pip
RUN python3 -m pip install -r requirements.txt

ENV port 8501
CMD bash -c "python3 -m streamlit run app.py --server.port=${port} --server.address=0.0.0.0"
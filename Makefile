install:
	python -m pip install -r requirements.txt
	python -m pip install -e .

demo:
	python run.py demo --data-dir data --config configs/demo.yaml --output-dir results/demo

test:
	python -m pytest -q

docker-build:
	docker build -t drosophila-pirna-germplasm .

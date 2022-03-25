VERSION=v1.5

build:
	docker build -t dbest/coverage:latest -f bam_coverage.dockerfile . 

publish:
	docker tag dbest/coverage:latest dbest/coverage:$(VERSION)
	docker push dbest/coverage:$(VERSION)

lint:
	pylint bam_coverage.py


format:
	black -v --line-length 1000 bam_coverage.py


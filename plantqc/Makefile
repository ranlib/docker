all:
	docker build -t dbest/plantqc:latest .

publish:
	docker tag dbest/plantqc:latest dbest/plantqc:v7.2
	docker push dbest/plantqc:v7.2

clean:
	rm *~

conda:
	conda env export > plantqc.yml
	conda env create -f plantqc.yml


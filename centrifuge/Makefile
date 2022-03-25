all:
	docker build -t dbest/centrifuge:latest .

publish:
	docker tag dbest/centrifuge:latest dbest/centrifuge:v1.0.4
	docker push dbest/centrifuge:v1.0.4

clean:
	rm *~


raconrounds:
	docker build -t dbest/raconrounds:latest -f racon.dockerfile .

raconrounds_publish:
	docker tag dbest/raconrounds:latest dbest/raconrounds:v2.2
	docker push dbest/raconrounds:v2.2

racon:
	docker build -t dbest/racon:v1.4.21 .

clean:
	rm *~

.PHONY: build test race bench lint fmt snapshot clean

build:
	go build -o dist/admix ./cmd/admix
	go build -o dist/server ./cmd/server

test:
	go test -count=1 ./...

race:
	go test -race -count=1 ./...

bench:
	go test -bench=. -benchmem -run='^$$' ./...

lint:
	golangci-lint run

fmt:
	gofmt -w .

snapshot:
	goreleaser build --snapshot --clean

clean:
	rm -rf dist

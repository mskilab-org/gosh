package main

import (
	"os"

	"github.com/mskilab-org/gosh/internal/cli"
)

var version = "v0.0.0-dev"

func main() {
	os.Exit(cli.Main(os.Args[1:], os.Stdout, os.Stderr, version))
}

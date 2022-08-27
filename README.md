# Bioseqdb

## Docker

Docker can be used to quickly spin up a PostgreSQL instance with the extension installed. To do so, run `docker run -d --name bioseqdb -e POSTGRES_PASSWORD=password -p 5432:5432 unneon/bioseqdb`. This will create and launch a instance in the background, to which you can connect using hostname `localhost`, port `5432`, username `postgres`, password `password` and database `postgres`. The data volume is persistent across reboots, and the container can be started again using `docker start bioseqdb`.

## Development

To work on the extension, install Postgres (`apt install postgresql postgresql-server-dev-all` on Ubuntu, `pacman -S postgresql` on Arch) and CMake, and set up a Postgres instance (see [Ubuntu docs](https://ubuntu.com/server/docs/databases-postgresql) or [Arch docs](https://wiki.archlinux.org/title/PostgreSQL)). To compile and install bioinformatics libraries used by the extension, follow the instructions from relevant steps in the [.github/workflows/ci.yml](.github/workflows/ci.yml) file.

To build and install the extension, create a `build/` directory and run `cmake ..` from it. You can now build the extension by running the `make` command in the build directory, and install it with `sudo make install`. A typical development flow is running `make && sudo make install && sudo systemctl restart postgresql`. The entire process should take about a second.

After first installing the extension, you need to run `CREATE EXTENSION bioseqdb;` to load the extension to the active database. If you modify the definitions of any SQL functions or types, remember to drop any affected tables, `DROP EXTENSION bioseqdb CASCADE;` and repeat the `CREATE EXTENSION` command.

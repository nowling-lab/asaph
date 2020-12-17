# Testing Asaph
All commits to the Asaph repository trigger Travis CI to run a series of integration tests.  If you are modifying the Asaph source code, it can be helpful to run the tests yourself to verify that your changes didn't break the code in any obvious ways.

## Running the Bats Tests
If you want to test the install, you can use the `bats` framework to run integration tests.  On Debian and Ubuntu systems, bats can be installed with apt:

```bash
$ sudo apt install bats
```

For other platforms, see the bats [documentation](https://bats-core.readthedocs.io/en/latest/).  Once installed, you can test your installation with:

```bash
$ bats bats_tests/*.bats
```

## What Next?
Now that Asaph is installed, check out some of our other [tutorials](README.md).

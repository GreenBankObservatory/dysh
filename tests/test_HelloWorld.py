class TestHelloWorld:
    value = 1

    def test_value(self):
        print("Hello World!")
        assert self.value == 1
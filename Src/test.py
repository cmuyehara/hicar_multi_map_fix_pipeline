import pysam

bam_file = "example.bam"


class MyRangeIterator:
    def __init__(self, start, end):
        self.current = start
        self.end = end

    def __iter__(self):
        # The __iter__ method should return the iterator object itself.
        return self

    def __next__(self):
        # The __next__ method returns the next item in the sequence.
        # It raises StopIteration when there are no more items.
        if self.current < self.end:
            value = self.current
            self.current += 1
            return value
        else:
            raise StopIteration

# Usage of the custom iterator
my_numbers = MyRangeIterator(1, 5)
while True:
    if any(my_numbers):
        number = next(my_numbers)
    else:
        break
    print(number)
    
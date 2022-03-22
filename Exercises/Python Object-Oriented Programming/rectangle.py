from math import sqrt


class Rectangle:
    def __init__(self, width, height):
        self.width = width
        self.height = height

    def area(self):
        return self.width * self.height

    def perimeter(self):
        return 2 * self.width + 2 * self.height

    def diagonal(self):
        return sqrt(self.width ** 2 + self.height ** 2)


class Square(Rectangle):
    def __init__(self, side):
        super().__init__(side, side)


if __name__ == "__main__":
    rectangles = [
        Rectangle(5, 10),
        Rectangle(5, 5),
        Square(5),
        Square(10)
    ]

    for i, rectangle in enumerate(rectangles):
        print(f"Rectangle {i + 1}")
        print(rectangle.area())
        print(rectangle.perimeter())
        print(rectangle.diagonal())
        print()

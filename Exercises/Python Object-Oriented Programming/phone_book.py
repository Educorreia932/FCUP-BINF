from copy import copy

class PhoneBook:
    def __init__(self, entries={}):
        self.entries = entries

    def add_entry(self, name, number):
        self.entries[name] = number

    def search_by_name(self, name):
        if name in self.entries:
            return self.entries[name]

    def search_by_number(self, number):
        for key, value in self.entries.items():
            if number == value:
                return key

    def __str__(self):
        return str(self.entries)

    def copy(self):
        return PhoneBook(copy(self.entries))

if __name__ == "__main__":
    phone_book = PhoneBook()

    phone_book.add_entry("Raul Aguilera", 919297989)
    phone_book.add_entry("Richard Thomas", 929296058)
    phone_book.add_entry("Quyen Mote", 938057559)
    phone_book.add_entry("Deborah Clark", 969324774)

    print(phone_book.search_by_name("Quyen Mote"))
    print(phone_book.search_by_number(969324774))
    print(phone_book)
    print(phone_book.copy())

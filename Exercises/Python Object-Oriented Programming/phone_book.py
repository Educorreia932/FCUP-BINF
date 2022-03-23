from copy import copy
from rich import print

class PhoneBook:
    def __init__(self, entries=[]):
        self.entries = entries

    def add_entry(self, name, number):
        self.entries.append((name, number))

    def search_by_name(self, name):
        for entry in self.entries:
            if entry[0] == name:
                return entry

    def search_by_number(self, number):
        for entry in self.entries:
            if entry[1] == number:
                return entry

    def __str__(self):
        return str(self.entries)

    def copy(self):
        return PhoneBook(copy(self.entries))

class EmailBook(PhoneBook):
    def __init__(self, entries=[]):
        super().__init__(entries)

    def add_entry(self, name, number, email):
        self.entries.append((name, number, email))

    def search_by_email(self, email):
        for entry in self.entries:
            if len(entry) > 2 and entry[2] == email:
                return entry


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

    email_book = EmailBook(phone_book.entries)

    email_book.add_entry("Paul Smith", 929123719, "paul@smith.com")

    print(email_book)

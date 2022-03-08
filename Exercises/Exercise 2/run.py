def statistics(l):
    l = list(filter(lambda x: x > 0, l))

    return {"min:": min(l), "max:": max(l), "avg:": sum(l) / len(l)} 

print(statistics([4, 5, -6, 6, 10, -15, 0, 5]))

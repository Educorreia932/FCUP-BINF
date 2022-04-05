## Task 1 – Regular Expressions

### 1.Describe a regular expression for the following patterns:

> AATATC, AAATATC, AAATATC  

```re
A{2,3}TATC
```

> ACGT, ACCT, ACTT, ACAT  

```re
AC[GCTA]T
```

> "She said heed very clearly!", "She Said head very clearly", "She said hood clearly"

```re
She [sS]aid h[ee][ea][oo]d (very)? clearly!?
```

### 2. Which of the following regular expressions are syntactically incorrect? Why? What kind of strings do the valid ones match?

> `?.\b`

Invalid

> `[?].\b`

Valid

> `[^\D]{4,1}`

Invalid

> `^[^^]`

Invalid

> `++`

Invalid

> `+\+`

Invalid

> `\++`

Valid

> `\+\+`

Valid

### 3. What kind of strings will be matched by the following regular expressions? In your answer, break them up into parts and explain what each one does.

> ^[+-]?\d+(\.\d+)?

Floating point and integer numbers

**Ex.:** 3.14, 21

> \b[aeiou][a-z]{,4}\b
> [.?!]\s+([A-Z][A-Za-z]*)

### 4. Write a regular expression that describes an animal name with three letters and the middle letter is a vowel.

```re
[a-z][aeiou][a-z]
```

### 6. Write an RE to capture the following sentences:

> Red. Blue. Green!
> Red. blue. Green!
> Red. Yellow. Blue. Green.
> Red. Blue. Grey.

```re
Red\.([bB]lue|Yellow)\.(Blue\.)?(Grey|Green)[!.] 
```

### 7.Which of the following statements are incorrect?

> a) ‘pq*’ will match ‘pq’

Yes

> b) ‘pq?’ matches ‘p’

Yes

> c) ‘p{4}, q’ does not match ‘pppq’

Yes

> d) ‘pq+’ matches ‘p’

No

### 8. Write a RE to match an IP address.

[0-255]\.[0-255]\.[0-255]\.[0-255]

### 9. Write a RE to match an email address.

```re
\w+@\w+\.\w+
```



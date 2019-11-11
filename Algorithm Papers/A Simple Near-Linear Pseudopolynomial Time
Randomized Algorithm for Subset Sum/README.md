#### Given a set contains ```m``` kinds of numbers, the number of ```a[i]``` is ```b[i]```, ```(b[i] != 0)```
#### Return there are number of such subsets S modulo ```p``` that subsum is ```[1. n]```, p is a prime

#### Example:
```
n = 5, m = 3
a[0] = 5  b[0] = 2
a[1] = 1  b[1] = 2
a[2] = 2  a[2] = 3
answer = 1, 2, 1, 2, 1 
```

#### explanation
```
1: {1}
2: {1, 1}, {2}
3: {1, 2}
4: {2, 2}, {1, 1, 2}
5: {1, 2, 2}
```
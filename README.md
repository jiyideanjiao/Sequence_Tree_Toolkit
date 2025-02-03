# Moving Ortholog

- Author: Chao Tong
- Date: Jan-1-2020


```
cat list | awk '{print$1}' | xargs -I '{}' mv {} yes/
cat list | awk '{print$1}' | xargs -I '{}' cp {} yes/
cat final | awk '{print$1}' | xargs -I '{}' mv {} yes/
```

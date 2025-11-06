# fastascan

 - rust enabled fastascan
 - algorithm that it uses for binning: How it classifies bins:
```
let max_len = *self.lengths.iter().max().unwrap_or(&0);
let bin_size = (max_len as f64 / 10.0).max(1.0) as usize;
let num_bins = (max_len + bin_size - 1) / bin_size + 1;
```

```
 cargo build
``` 

Gaurav Sablok \
codeprog@icloud.come

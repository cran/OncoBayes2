  /*
   * Author: S. Weber
   *
   * Various utlities
   */
  
  // turn a slicing variable for a ragged array
  // S = {5, 6, 11}
  // into
  // Si = {0, 5, 5 + 6, 5+6+11} + 1
  // such that we can index the ragged array A as
  // A[Si[u] : Si[u+1]-1]
  // for the uth unit
  array[] int make_slice_index(array[] int S) {
    array[size(S) + 1] int Si;
    int cv = 1;
    Si[1] = cv;
    for (i in 1 : size(S)) {
      cv = cv + S[i];
      Si[i + 1] = cv;
    }
    return Si;
  }
  
  // helper function for rle_int
  int rle_elem_count(array[] int set) {
    int U = 1;
    for (i in 2 : num_elements(set)) {
      if (set[i - 1] != set[i]) {
        U = U + 1;
      }
    }
    return U;
  }
  
  // helper function for rle_int
  int rle_elem_count_vector(vector set) {
    int U = 1;
    for (i in 2 : num_elements(set)) {
      if (set[i - 1] != set[i]) {
        U = U + 1;
      }
    }
    return U;
  }
  
  // repeated length encoding, see rle in R
  // this function returns the number of times elements are repeated ...
  array[] int rle_int(array[] int set) {
    array[rle_elem_count(set)] int res;
    int c = 1;
    res[1] = 1;
    for (i in 2 : num_elements(set)) {
      if (set[i - 1] == set[i]) {
        res[c] = res[c] + 1;
      } else {
        c = c + 1;
        res[c] = 1;
      }
    }
    return res;
  }
  
  // ... and this function returns which elements are repeated for ints
  array[] int rle_elem_int(array[] int set) {
    int N = rle_elem_count(set);
    array[N] int first_ind = make_slice_index(rle_int(set))[1 : N];
    
    return set[first_ind];
  }
  
  array[] int decimal2base(int decimal, int digits, int base) {
    array[digits] int base_rep;
    int current = decimal;
    for (i in 1 : digits) {
      base_rep[i] = current % base;
      current = current %/% base;
    }
    return base_rep;
  }
  
  int power_int(int number, int power) {
    if (power < 0) {
      reject("Cannot raise an integer to a negative power and expect an integer result.");
    }
    if (power == 0) {
      return 1;
    } else {
      return number * power_int(number, power - 1);
    }
  }
  
  // count the number of unique elements
  int cardinality_int(array[] int elems) {
    return rle_elem_count(sort_asc(elems));
  }
  
  // count the number of unique elements
  int cardinality_vector(vector elems) {
    vector[num_elements(elems)] sort_asc_elems = sort_asc(elems);
    return rle_elem_count_vector(sort_asc_elems);
  }
  
  // create an integer sequence
  array[] int seq_int(int start, int end) {
    int N = end - start + 1;
    array[N] int seq;
    for (i in 1 : N) {
      seq[i] = i + start - 1;
    }
    return seq;
  }
  
  // return an int vector with each element in the input repeated as
  // many times as given in the each argument
  array[] int rep_each(array[] int set, array[] int each) {
    int N = sum(each);
    array[N] int replicated;
    int p = 1;
    
    for (i in 1 : size(set)) {
      replicated[p : p + each[i] - 1] = rep_array(set[i], each[i]);
      p = p + each[i];
    }
    
    return replicated;
  }
  
  
  // count number times elem appears in test set
  int count_elem(array[] int test, int elem) {
    int count;
    count = 0;
    for (i in 1 : num_elements(test)) {
      if (test[i] == elem) {
        count = count + 1;
      }
    }
    return count;
  }
  
  // count number times elems appears in test set
  array[] int count_elems(array[] int test, array[] int elems) {
    array[num_elements(elems)] int counts;
    for (i in 1 : num_elements(elems)) {
      counts[i] = count_elem(test, elems[i]);
    }
    return counts;
  }
  
  // find elements in test which are equal to elem
  array[] int which_elem(array[] int test, int elem) {
    array[count_elem(test, elem)] int res;
    int ci;
    ci = 1;
    for (i in 1 : num_elements(test)) {
      if (test[i] == elem) {
        res[ci] = i;
        ci = ci + 1;
      }
    }
    return res;
  }
  
  

#!/usr/bin/env python3

import numpy as np
from itertools import chain, combinations

class test_yield_func():
    def __init__(self, seed, iter_val):
        self.seed = seed
        self.iter_val = iter_val

    def powerset_generator(self, num_passing):
        '''create an iterator object which makes a powerset, excluding the first empty set and the last full set
        By definition, powerset([1,2,3]) --> () (1,) (2,) (3,) (1,2) (1,3) (2,3) (1,2,3)
        In practice, powerset_generator([1,2,3]) --> (1,) (2,) (3,) (1,2) (1,3) (2,3)
        Further, because it's a generator and I'll use next() to grab subsets (while I want to randomly pick subsets) I've got some fancy shuffling.
        Sadly, that means all of the subsets of a certain length will not be shuffled with respect to each other'''
        np.random.seed(num_passing//10)
        s = np.random.permutation(num_passing) #a shuffled np.arange(num_passing), these are the elements that will be combined
        t = np.arange(1, num_passing) #these are the elements that will be used for the length of the combinations, we start with 1 to avoid the empty set
        np.random.shuffle(t) #we randomly shuffle t so that we don't always get all the 1 element subsets and then all the 2 element subsets etc.
        powerset_iterator = chain.from_iterable(combinations(s, r) for r in t)
        #print(list(powerset_iterator))
        return powerset_iterator


    def find_num_in_powerset(self, num_passing):
        return ((2**num_passing)-2)

    def move_iterator_for_sorted(self, iterator_obj, current_i, next_i):
        while True:
            try:
                new_subset = next(iterator_obj)
                current_i += 1
            except StopIteration:
                '''shouldn't actually need this'''
                logging.error("Something bad happened with iteration")
                break
            if current_i == next_i:
                new_subset = next(iterator_obj)
                break
            elif current_i > next_i:
                break
        yield new_subset

    def indice_iterator(self, powerset_size, gene_index):
        np.random.seed(int((gene_index+1)*self.seed))
        indice_iterator = iter(sorted(np.random.choice(np.arange(powerset_size), size=self.iter_val, replace=False)))
        #print(list(indice_iterator))
        return indice_iterator

    def drive_random_iter(self):
        num_passing =  4
        gene_index = 1
        powerset_it_obj = self.powerset_generator(num_passing)
        powerset_size = self.find_num_in_powerset(num_passing)
        indice_it_obj = self.indice_iterator(powerset_size, gene_index)
        current_i = 0
        while True:
            print("CI: ", current_i)
            try:
                next_indice = next(indice_it_obj)
                print("NI: ", next_indice)
            except StopIteration:
                break
            new_subset = np.array(list(self.move_iterator_for_sorted(powerset_it_obj, current_i, next_indice)))
            current_i = next_indice + 1
            print("NS: ", new_subset)


test_yield_instance = test_yield_func(42, 5)
test_yield_instance.drive_random_iter()

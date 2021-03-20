from collections import Counter


labels = []
sorted_list = []
occurrence = []



class Watershed:
    class Index:
        """
        purpose: object to hold start and stop indices of neighborhood
        """

        def __init__(self, start, stop):
            self.start = start
            self.stop = stop

    def __init__(self, data, threshold):
        """
        :param data: list containing expression data
        """

        self.data = data
        self.threshold = threshold
        self.abs_data = [abs(val) for val in self.data]  # to find positive & negative peaks

    def init_global(self):
        """
        :return:
        1. labels - List to peak label of each data point
        2. sorted_list - absolute expression values sorted in reverse
        """
        global labels
        global sorted_list

        labels = [0] * len(self.data)
        sorted_list.extend(self.abs_data)
        sorted_list.sort(reverse=True)

    def current_max(self):
        """
        :return: maximum from unlabelled data points
        """

        global sorted_list
        global occurrence

        val = sorted_list[0]

        if self.abs_data.count(val) == 1:
            sorted_list.pop(0)
            return self.abs_data.index(val)

        else:  # if val occurs multiple times

            if not occurrence:
                occurrence = [i for i, x in enumerate(self.abs_data) if
                              x == val]  # indices of occurrences of current val

            idx = occurrence[0]
            sorted_list.pop(0)
            occurrence.pop(0)

            return idx

    def adj_range(self, index, threshold):
        """
        :param index: index of current datapoint
        :param threshold: neighborhood threshold
        :return: start and stop indices of adjacency list (list of neighbors within given threshold)
        """

        i = 0
        j = 0
        start = index
        stop = index
        while i < threshold and start > 0:  # ensure start doesnt go below the 0th index
            start -= 1
            i += 1

        while j < threshold and stop < len(self.data):  # ensure stop doesnt go above list size
            stop += 1
            j += 1

        return Watershed.Index(start, stop)

    def adj_list(self, index, start, stop):
        """
        :param index: index of current datapoint
        :param start: start index of adjacency list
        :param  stop: stop index of adjacency list
        :return: list of neighbors of current datapoint
        """

        adj_list = self.abs_data[start:stop + 1]
        adj_list.pop(index - start)  # remove current datapoint from neighbors list
        return adj_list

    def closest_neighbor(self, index, start, stop, largest_neighbor):
        """
        :param index: index of current datapoint
        :param start: start index of adjacency list
        :param stop: stop index of adjacency list
        :param largest_neighbor: largest datapoint in neighborhood (adjacency list)
        :return: index of largest neighbor occurring closest to current data point
        """

        # get indices of all occurrences of largest neighbor from neighborhood[start:stop)
        neighbor_idx = [i for i, x in enumerate(self.abs_data[start:stop + 1]) if x == largest_neighbor]

        for i in range(len(neighbor_idx)):  # adjust indices with offset for the whole list
            neighbor_idx[i] += start

        neighbor_idx = [abs(x - index) for x in neighbor_idx]  # distance from current datapoint

        return neighbor_idx.index(min(neighbor_idx))

    def closest_labelled_neighbor(self, index, start, stop, largest_neighbor):
        """
        :param index: index of current datapoint
        :param start: start index of adjacency list
        :param stop: stop index of adjacency list
        :param largest_neighbor: largest datapoint in neighborhood (adjacency list)
        :return: index of largest labelled neighbor occurring closest to current data point
        """
        global labels

        # get indices of all occurrences of largest neighbor from neighborhood[start:stop)
        neighbor_idx = [i for i, x in enumerate(self.abs_data[start:stop + 1]) if x == largest_neighbor]

        for i in range(len(neighbor_idx)):  # adjust indices with offset for the whole list
            neighbor_idx[i] += start

        labelled_idx = []
        for i in range(len(neighbor_idx)):  # find indices of labelled neighbors
            if labels[neighbor_idx[i]] != 0:
                labelled_idx.append(neighbor_idx[i])

        if len(labelled_idx) != 0:
            labelled_dist = [abs(x - index) for x in labelled_idx]  # closest labelled neighbor to current datapoint

            return labelled_idx[labelled_dist.index(min(labelled_dist))]  # return closest labelled neighbor

        else:  # if no neighbors are labelled

            return -1

    def watershed(self):
        """
        #specifies which distance between points is considered near enough
        :return: every datapoint is labelled with appropriate peak identifier
        """

        global labels

        i = 0  # keep track of labelled data points
        peak = 1
        while i < len(self.abs_data):

            global labels
            local_labels = []
            local_labels.extend(labels)

            max_idx = self.current_max()  # find the index of largest unlabelled data point

            neighborhood = self.adj_range(max_idx, self.threshold)  # get indices of adjacency list
            adj_list = self.adj_list(max_idx, neighborhood.start, neighborhood.stop)
            largest_neighbor = max(adj_list)

            if largest_neighbor != self.abs_data[max_idx]:  # if largest neighbor is not equal to current value

                if adj_list.count(largest_neighbor) == 1:  # single occurrence of largest neighbor

                    if largest_neighbor > self.abs_data[max_idx]:  # give same label as largest neighbor

                        neighbor_idx = self.abs_data.index(largest_neighbor, neighborhood.start,
                                                           neighborhood.stop + 1)  # index of largest neighbor in
                        # adjacency list

                        # for item in range(0, max_idx):
                        labels[max_idx] = labels[neighbor_idx]

                    else:  # give new label

                        # for item in range(0, max_idx + 1):
                        labels[max_idx] = peak
                        peak += 1

                else:  # multiple occurrences of largest neighbor

                    if largest_neighbor > self.abs_data[max_idx]:  # give same label as largest neighbor

                        neighbor_idx = self.closest_neighbor(max_idx, neighborhood.start,
                                                             neighborhood.stop, largest_neighbor)

                        labels[max_idx] = labels[neighbor_idx]

                    else:  # give new label

                        labels[max_idx] = peak
                        peak += 1

            else:  # if largest neighbor is equal to current value

                if adj_list.count(largest_neighbor) == 1:  # single occurrence of largest neighbor

                    neighbor_idx = self.abs_data.index(largest_neighbor, neighborhood.start,
                                                       neighborhood.stop)  # index of largest neighbor in adjacency list
                    labels[max_idx] = labels[neighbor_idx]

                else:  # multiple occurrences of largest neighbor

                    neighbor_idx = self.closest_labelled_neighbor(max_idx, neighborhood.start,
                                                                  neighborhood.stop, largest_neighbor)

                    if neighbor_idx != -1:  # if neighbor is labelled

                        labels[max_idx] = labels[neighbor_idx]

                    else:  # give new label

                        labels[max_idx] = peak
                        peak += 1

            i += 1

    def peak_calling(self):

        peaks = []
        neighborhoods = [k for k, cnt in Counter(labels).items() if
                         cnt > 1]  # get unique identifier of each neighborhood

        for neighborhood_label in neighborhoods:
            current_neighborhood_idx = []
            for i, x in enumerate(labels):
                if x == neighborhood_label:
                    current_neighborhood_idx.append(i)  # indices of current neighborhood

            peaks.append(max(self.data[current_neighborhood_idx[0]: current_neighborhood_idx[-1] + 1],
                             key=abs))  # peak of current neighborhood

        return peaks



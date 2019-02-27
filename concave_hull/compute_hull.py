from scipy.spatial import ConvexHull
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull

class ConcHuller(object):

    def __init__(self, n_param, points, debug = False, savedir = './'):

        """
        assume that everything is in rectangular coordinate axis
        FIXME: add axis converter later on
        FIXME: points - row: data, col: x, y, z, ...
        """

        self.n_param = n_param
        self.points = points
        self.debug = debug
        self.savedir = savedir

        if self.points.shape[1] == 2: # two dimensions
            self.two_dim = True
        else:
            self.two_dim = False

    def _get_vec_mag(self, vec):
        return np.sqrt(np.dot(vec, vec))

    def _get_simplex_param(self, pts):
        if self.two_dim:
            vec = pts[1, :] - pts[0, :]
            #print ("here now:", pts, self._get_vec_mag(vec))
            return (self._get_vec_mag(vec))
        else:
            mag0 = self._get_vec_mag(pts[1, :] - pts[0, :])
            mag1 = self._get_vec_mag(pts[2, :] - pts[0, :])
            mag2 = self._get_vec_mag(pts[2, :] - pts[1, :])

            return((mag0 + mag1 + mag2)/3.)

    def _get_min_dist_pt_simp_pts(self, pt, simplex_pts):

        vec = pt - simplex_pts
        return (np.min(np.sqrt(np.diag(np.matmul(vec, vec.T))))) # bad readability??


    # should not use normal in computing for the distance as far points highly acute to the edge would be closer
    def _get_normal(self, pts):

        # 2D case
        if self.two_dim:
            vec = pts[1, :] - pts[0, :]
            norm = np.array([-1*vec[1], vec[0]])
            return(norm/self._get_vec_mag(norm)) # more readable when using rotation matrix?
        # 3D case
        else:
            vec0 = pts[1, :] - pts[0, :]
            vec1 = pts[2, :] - pts[0, :]
            norm = np.cross(vec0, vec1)
            return(norm/self._get_vec_mag(norm)) # no need to worry direction of normal vector

    def compute_hull(self):

        # compute convex hull
        self.temp_hull = ConvexHull(self.points)

        self.hull_vertices = self.temp_hull.vertices # update

        # prioritize breaking up of long segments
        temp_param = np.array([self._get_simplex_param(self.points[x]) for x in self.temp_hull.simplices])
        sort_idx = np.argsort(temp_param) # index for smallest simplex param to largest
        sort_idx = sort_idx[::-1] # reverse order (largest to smallest) -> prioritize large parameter values
        self.hull_simplices = self.temp_hull.simplices[sort_idx] # rearrange (largest to smallest parameter value)
        #print(temp_param)

        non_vert = np.array([x for x in range(self.points.shape[0])])
        mask = np.ones(self.points.shape[0], dtype=bool)
        mask[self.hull_vertices] = False
        non_vert = non_vert[mask]

        idx = 0
        rm_simplices = []
        while True:
            print("iteration: ", idx)
            #for x in range(self.hull_simplices.shape[0]):
            #    print(x, self.hull_simplices[x])
            #print("end of simplex list")
            simplex_idx = self.hull_simplices[idx]
            #print("simplex index: ", simplex_idx)
            # norm = self._get_normal(self.points[simplex_idx]) # using norm is wrong!!!
            simp_centroid = np.mean(self.points[simplex_idx], axis = 0) # average along rows
            #print("simplex centroid: ", simp_centroid)

            # compute distance between current simplex centroid and non-edges - or should you use nearest simplex point?
            sep = simp_centroid - self.points[non_vert]
            dist = np.sqrt(np.diag(np.matmul(sep, sep.T)))

            point_idx = np.argsort(dist) # idx of nearest point to simplex from nearest to farthest

            #print("candidate point_idx: ", self.points[non_vert[point_idx[0]]])
            # check neighbor nodes only once - if other edges are also near to the point, then reject it
            found = False
            nearest_idx = point_idx[0]
            #print("nearest: ", nearest_idx)
            centroids = np.mean(self.points[self.hull_simplices], axis = 1)
            sep = centroids - self.points[non_vert[point_idx[0]]]
            dist2 = np.sqrt(np.diag(np.matmul(sep, sep.T)))
            edge_idx = np.argsort(dist2)

            if edge_idx[0] == idx: # if edge nearest to point is the current edge (idx)
                nearest_pt = self.points[non_vert[nearest_idx]]
                #print("nearest pt: ", nearest_pt)
                found = True
            ##### END OF ONE-TIME CHECK#######

            """
            # for ITERATIVELY checking neighbor edges - results to too much digging
            found = False
            for  x in point_idx:
                #print("x", x)
                nearest_idx = non_vert[x] # candidate point nearest to edge
                #print("nearest_idx: ", nearest_idx)
                #print("hull simplices:", self.hull_simplices)
                centroids = np.mean(self.points[self.hull_simplices], axis = 1)
                #print("centroids: ", centroids)
                #print("shape: ", centroids.shape)
                #print("shape2: ", self.points[self.hull_simplices].shape)
                sep = centroids - self.points[nearest_idx]
                dist2 = np.sqrt(np.diag(np.matmul(sep, sep.T)))
                edge_idx = np.argsort(dist2)

                #print("here", dist2)
                if edge_idx[0] == idx: # if edge nearest to point is the current edge (idx)
                    nearest_pt = self.points[non_vert[x]]
                    found = True
                    break # found the point nearest to the edge
            """

            # if did not find anything that is uniquely near to current edge
            if not found:
                #print("found nothing")
                idx += 1
                if (idx == self.hull_simplices.shape[0]):
                    break # already done
                else:
                    continue # skip loop

            simplex_param = self._get_simplex_param(self.points[simplex_idx])

            sep = nearest_pt - self.points[simplex_idx]
            #print("separation matrix: ", sep)
            pt_simp_sep = (np.min(np.sqrt(np.diag(np.matmul(sep, sep.T))))) # min distance from point to simplex pts
            #print("pt to simplex separation: ", pt_simp_sep)
            #print("ratio: ", pt_simp_sep/simplex_param)

            #if (simplex_param/pt_simp_sep) > self.n_param:
            #print(pt_simp_sep, simplex_param)
            if (pt_simp_sep/simplex_param) <= self.n_param:
                if self.two_dim:
                    to_add = np.array([[simplex_idx[0], non_vert[nearest_idx]], [simplex_idx[1], non_vert[nearest_idx]]])
                else:
                    to_add = np.array([[simplex_idx[0], simplex_idx[1], non_vert[nearest_idx]], [simplex_idx[1], simplex_idx[2], non_vert[nearest_idx]], [simplex_idx[2], simplex_idx[0], non_vert[nearest_idx]]])

                #print("appending: ", nearest_idx)
                self.hull_vertices = np.append(self.hull_vertices, non_vert[nearest_idx])

                self.hull_simplices = np.vstack((self.hull_simplices, to_add))

                rm_simplices.append(idx) # index of simplex to remove
                #print("before mask: ", non_vert)
                mask = (non_vert != non_vert[nearest_idx]) # remove point that is now a vertex of a simplex
                non_vert = non_vert[mask]
                #print("after mask: ", non_vert)

                # sort hull_simplices in idx + 1 to end - cool move
                temp_simplex = self.hull_simplices[idx + 1:]
                temp_param = np.array([self._get_simplex_param(self.points[x]) for x in temp_simplex])
                sort_idx = np.argsort(temp_param) # index for smallest simplex param to largest
                sort_idx = sort_idx[::-1] # reverse order (largest to smallest) -> prioritize large parameter values
                self.hull_simplices = np.vstack((self.hull_simplices[: idx + 1], temp_simplex[sort_idx])) # rearrange (largest to smallest parameter value)

                if self.debug and self.two_dim:
                    plt.plot(self.points[:, 0], self.points[:, 1], 'o')
                    for i in range(self.hull_simplices.shape[0]):
                        if i not in rm_simplices:
                            plt.plot(self.points[self.hull_simplices[i], 0], self.points[self.hull_simplices[i], 1], 'k-')

                    for x, y in to_add:
                        #plt.plot(self.points[x, 0], self.points[x, 1], 'r-', linewidth = 3)
                        #plt.plot(self.points[y, 0], self.points[y, 1], 'r-', linewidth = 3)
                        #print("annotate: ", self.points[x], self.points[y])
                        plt.annotate('O', self.points[x], fontsize = 12)
                        plt.annotate('O', self.points[y], fontsize = 12)


                    plt.savefig(self.savedir + '/' + str(idx) + '.jpg')
                    plt.clf()
                    #print()

            idx += 1
            if (idx == self.hull_simplices.shape[0]) or len(non_vert) == 0:
                break

        mask = np.ones(self.hull_simplices.shape[0], dtype=bool)
        mask[rm_simplices] = False
        self.hull_simplices = self.hull_simplices[mask]

        return

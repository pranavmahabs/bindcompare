# Python Code for Interval tree
class BindingSite:
    def __init__(self, low, high):
        self.low = low
        self.high = high
 
    def __str__(self):
        return "[" + str(self.low) + "," + str(self.high) + "]"
 
 
class Node:
    def __init__(self, range, max):
        self.range = range
        self.max = max
        self.left = None
        self.right = None
        # self.protein = name
 
    def __str__(self):
        return "[" + str(self.range.low) + ", " + str(self.range.high) + "] " + "max = " + str(self.max) + "\n"
 
 
def insert(root, x):
    if root == None:
        return Node(x, x.high)
 
    if x.low < root.range.low:
        root.left = insert(root.left, x)
    else:
        root.right = insert(root.right, x)
 
    if root.max < x.high:
        root.max = x.high
 
    return root
 
 
def inOrder(root):
    if root == None:
        return
 
    inOrder(root.left)
    print(root, end="")
    inOrder(root.right)

def search(root, x):
    if root == None:
        # return a dummy interval range
        return []
    overlaps = []
 
    # if x overlaps with root's interval
    rlow = root.range.low
    rhigh = root.range.high
    print(root)
    if (x.low > root.max):
        return overlaps
    elif (rlow <= x.high and rhigh >= x.low):
        overlaps.append(root.range.__str__())
 
    # Appropriately Recur through the tree
    left_exists = root.left != None 
    if(left_exists and x.low < rlow and x.high < rlow):
        if (x.low > root.left.max):
            return overlaps
        else:
            return overlaps + search(root.left, x)
    elif(root.right != None and left_exists and x.low > root.left.max):
        if (x.low > root.right.max):
            return overlaps
        else:
            return overlaps + search(root.right, x)
    else:
        if (left_exists and root.left.max > x.low):
            overlaps.extend(search(root.left, x))
        overlaps.extend(search(root.right, x))
        return overlaps
 
 
if __name__ == '__main__':
    root = None
    root = insert(None, BindingSite(15, 20))
    root = insert(root, BindingSite(10, 30))
    root = insert(root, BindingSite(17, 19))
    root = insert(root, BindingSite(5, 20))
    root = insert(root, BindingSite(12, 15))
    root = insert(root, BindingSite(30, 40))
 
    print("Inorder traversal of constructed Interval Tree is")
    inOrder(root)
    print()
    i = BindingSite(15, 20)
    print("Searching for interval", i)
    print("Overlaps with ", search(root, i))
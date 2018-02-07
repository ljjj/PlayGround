/**
 * Definition for a binary tree node.
 * struct TreeNode {
 *     int val;
 *     TreeNode *left;
 *     TreeNode *right;
 *     TreeNode(int x) : val(x), left(NULL), right(NULL) {}
 * };
 */
class Solution {
private:
    bool isLeaf(TreeNode* root) {
        if (root == NULL) return false;
        return root->left == NULL && root->right == NULL;
    }
    
    bool houseRobbed(TreeNode* root) {
        if (root == NULL) return false;
        if (isLeaf(root)) return true;
        int left = 0;
        if (root->left != NULL) left = root->left->val;
        int right = 0;
        if (root->right != NULL) right = root->right->val;
        return root->val != left + right;
    }
    
public:
    int rob(TreeNode* root) {
        if (root == NULL) return 0;
        if (isLeaf(root)) return root->val;
        
        if (root->left == NULL && isLeaf(root->right)) {
            root->val = max(root->val, root->right->val);
            return root->val;
        }
        
        if (root->right == NULL && isLeaf(root->left)) {
            root->val = max(root->val, root->left->val);
            return root->val;
        }
        
        if (isLeaf(root->left) && isLeaf(root->right)) {
            root->val = max(root->val, root->left->val + root->right->val);
            return root->val;
        }
        
        if (root->left == NULL) {
            int robRight = rob(root->right);
            int rightLeft = 0;
            if (root->right->left != NULL) rightLeft = root->right->left->val;
            int rightRight = 0;
            if (root->right->right != NULL) rightRight = root->right->right->val;
            
            int robRoot = root->val + rightLeft + rightRight;
            if (houseRobbed(root->right)) root->val = max(robRoot, root->right->val);
            else root->val = robRoot;
            return root->val;
        }
        
        if (root->right == NULL) {
            int robLeft = rob(root->left);
            int leftLeft = 0;
            if (root->left->left != NULL) leftLeft = root->left->left->val;
            int leftRight = 0;
            if (root->left->right != NULL) leftRight = root->left->right->val;
            
            int robRoot = root->val + leftLeft + leftRight;
            if (houseRobbed(root->left)) root->val = max(robRoot, root->left->val);
            else root->val = robRoot;
            return root->val;
        }
        
        int robLeft = rob(root->left);
        int leftLeft = 0;
        if (root->left->left != NULL) leftLeft = root->left->left->val;
        int leftRight = 0;
        if (root->left->right != NULL) leftRight = root->left->right->val;
        
        int robRight = rob(root->right);
        int rightLeft = 0;
        if (root->right->left != NULL) rightLeft = root->right->left->val;
        int rightRight = 0;
        if (root->right->right != NULL) rightRight = root->right->right->val;
        
        int robRoot = root->val + leftLeft + leftRight + rightLeft + rightRight;
        if (houseRobbed(root->left) || houseRobbed(root->right)) root->val = max(robRoot, root->left->val + root->right->val);
        else root->val = robRoot;
        return root->val;
    }
};
# def create_dataset(data, n_predictions, n_next):
#     '''
#     对数据进行处理
#     '''
#     dim = data.shape[1]
#     train_X, train_Y = [], []
#     for i in range(data.shape[0] - n_predictions - n_next):
#         # for i in range(data.shape[0] - n_predictions - n_next - 1):
#         if i==1:
#             a = data[i:(i + n_predictions), :]
#             train_X.append(a)
#         else:
#             a = data[i:(i + n_predictions - 3), :]
#             train_X.append(a)
#         tempb = data[(i + n_predictions):(i + n_predictions + n_next), :]
#         b = []
#         for j in range(len(tempb)):
#             for k in range(dim):
#                 b.append(tempb[j, k])
#         train_Y.append(b)
#     train_X = np.array(train_X, dtype='float64')
#     train_Y = np.array(train_Y, dtype='float64')


x, y = [], []
x += ['a']
print(x)
import numpy as np
from sklearn.cross_decomposition import PLSRegression
from sklearn.svm import SVC
from collections import defaultdict

class CompositeModel:
  def __init__(self, gexp_matrix, major_type, subtype, gene_names):
    self.major_types = []
    self.subtypes = {}
    self.gene_names = gene_names
    exclude = ["unknown", "NA"]
    for i in range(len(major_type)):
      if major_type[i] not in self.major_types:
        self.major_types.append(major_type[i])
        self.subtypes[major_type[i]] = []
      if subtype[i] not in exclude and subtype[i] not in self.subtypes[major_type[i]]:
        self.subtypes[major_type[i]].append(subtype[i])
    self.subtype_list = sorted(list(set([s for t in self.major_types for s in self.subtypes[t]])))

    # build submodels to go into composite model
    self.major_submodels = []
    self.subtype_submodels = defaultdict(list)
    # major type vs all else (not OTHER)
    for t in self.major_types:
      truth = np.array([(0 if m == t else 1) for m in major_type])
      for nc in range(5,11):
        model = PLSRegression(n_components=nc)
        model.fit(gexp_matrix, truth)
        self.major_submodels.append(model)
      # subtype vs all else
      for s in self.subtypes[t]:
        idx = np.array([l for l in range(len(subtype)) if major_type[l] == t and subtype[l] not in exclude])
        truth = np.array([(0 if s == subtype[l] else 1) for l in range(len(subtype)) if l in idx])
        for nc in range(5,11):
          model = PLSRegression(n_components=nc)
          model.fit(gexp_matrix[idx,:], truth)
          self.subtype_submodels[t].append(model)
      # subtype pairwise
      for i in range(len(self.subtypes[t])-1):
        for j in range(i+1, len(self.subtypes[t])):
          idx = np.array([l for l in range(len(subtype)) if major_type[l] == t and (self.subtypes[t][i] in subtype[l] or self.subtypes[t][j] in subtype[l])])
          truth = np.array([(0 if self.subtypes[t][i] in subtype[l] else 1) for l in idx])
          if truth.sum() == 0 or truth.sum() == truth.shape[0]: # all 0 or all 1
            continue
          for nc in range(5,11):
            model = PLSRegression(n_components=nc)
            model.fit(gexp_matrix[idx,:], truth)
            self.subtype_submodels[t].append(model)
    # major pairwise
    for i in range(len(self.major_types)-1):
      for j in range(i+1, len(self.major_types)):
        idx = np.array([l for l in range(len(major_type)) if major_type[l] == self.major_types[i] or major_type[l] == self.major_types[j]])
        truth = np.array([(0 if major_type[l] == self.major_types[i] else 1) for l in idx])
        if truth.sum() == 0 or truth.sum() == truth.shape[0]: # all 0 or all 1
          continue
        for nc in range(5,11):
          model = PLSRegression(n_components=nc)
          model.fit(gexp_matrix[idx,:], truth)
          self.major_submodels.append(model)

    # train composite random forest on output from PLSDA submodels
    major_composite_input = self.make_plsda_composite(gexp_matrix, subs=False)
    sub_composite_input = self.make_plsda_composite(gexp_matrix, subs=True)
    major_idx = [i for i in range(len(major_type)) if major_type[i] in self.major_types]
    major_truth = np.array([self.major_types.index(major_type[i]) for i in range(len(major_type)) if i in major_idx])
    self.composite_major = SVC(probability=True)
    self.composite_major.fit(sub_composite_input[major_idx,:], major_truth)
    self.composite_subtype = {}
    self.subtype_names = {}
    for t in self.major_types:
      # this check at the end makes sure the subtype is in the list with fuzzy matching vvvvv
      subtype_idx = [i for i in range(len(subtype)) if major_type[i] == t and len([s for s in self.subtype_list if s in subtype[i]]) > 0]
      if len(subtype_idx) == 0:
        continue
      self.subtype_names[t] = [self.subtype_list[j] for j in range(len(self.subtype_list)) if len([i for i in subtype_idx if self.subtype_list[j] in subtype[i]]) > 0]
      subtype_truth = np.array([[j for j in range(len(self.subtype_list)) if self.subtype_list[j] in subtype[i]][0] for i in range(len(subtype)) if i in subtype_idx])
      if len([a for a in subtype_truth if a == subtype_truth[0]]) == len(subtype_truth): # has no variation, SVM won't run, others won't matter
        continue
      self.composite_subtype[t] = SVC(probability=True)
      self.composite_subtype[t].fit(sub_composite_input[subtype_idx,:], subtype_truth)

  def make_plsda_composite(self, gexp_matrix, subs=True):
    composite_input = []
    for m in self.major_submodels:
      o = m.predict(gexp_matrix).flatten()
      composite_input.append(o)
    if subs:
      for t in self.major_types:
        for m in self.subtype_submodels[t]:
          o = m.predict(gexp_matrix).flatten()
          composite_input.append(o)
    return np.array(composite_input).transpose()

  def predict(self, gexp_matrix):
    major_composite_input = self.make_plsda_composite(gexp_matrix, subs=False)
    sub_composite_input = self.make_plsda_composite(gexp_matrix, subs=True)
    major_probs = self.composite_major.predict_proba(sub_composite_input)
    subtype_probs ={}
    for t in self.composite_subtype:
      subtype_probs[t] = self.composite_subtype[t].predict_proba(sub_composite_input)
    ret = []
    for i in range(major_probs.shape[0]):
      major_pred = self.major_types[np.argmax(major_probs[i])]
      ret.append(
        (
          major_pred,
          self.subtype_names[major_pred][np.argmax(subtype_probs[major_pred])] if major_pred in self.composite_subtype else None,
          major_probs,
          [subtype_probs[s] for s in subtype_probs]
        )
      )
    return ret

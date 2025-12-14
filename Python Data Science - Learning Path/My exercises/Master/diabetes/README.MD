Diabetes Prediction and Model Interpretability
Project Overview

This project focuses on predicting diabetes using supervised machine learning models, with an emphasis on model evaluation, comparison, and interpretability. The goal is not only to achieve strong predictive performance but also to understand how clinical features influence model decisions.

Objectives

Train and evaluate multiple machine learning models for diabetes prediction

Compare model performance using appropriate metrics for imbalanced classification

Identify the most predictive clinical features

Apply interpretable machine learning techniques to explain model behavior

Dataset and Preprocessing

The dataset contains clinical features relevant to diabetes diagnosis

Exploratory Data Analysis (EDA) was performed to understand feature distributions and correlations

The data was split into training and testing sets

Standard preprocessing steps were applied where necessary

Models Implemented

Logistic Regression (baseline model)

Random Forest

Gradient Boosting

Hyperparameter tuning was performed using RandomizedSearchCV for tree-based models.

Model Evaluation

Models were evaluated using multiple metrics, with special attention to class imbalance:

Confusion Matrix

Precision, Recall, and F1-score

ROC-AUC

Precision–Recall Curve

PR-AUC (Average Precision Score)

Model comparison showed that Gradient Boosting achieved the best overall performance in terms of PR-AUC and precision–recall trade-offs.

Model Interpretability

To ensure interpretability and clinical relevance, the following techniques were applied:

Permutation feature importance to identify key predictors

Partial Dependence Plots (PDP) to visualize the effect of important features on predicted probability

In particular, the Partial Dependence Plot for HbA1c level demonstrated a sharp increase in predicted diabetes probability beyond clinically established thresholds, aligning with medical knowledge.

Conclusions

Gradient Boosting was selected as the final model based on superior PR-AUC performance

HbA1c level emerged as one of the most influential predictors

Model behavior is consistent with established clinical definitions of diabetes

Interpretability techniques provided valuable insights into feature effects

Future Work

Explore additional boosting models such as XGBoost or LightGBM

Apply SHAP values for more granular interpretability

Evaluate model performance across different clinical subgroups
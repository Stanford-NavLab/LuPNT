import argparse

import mlflow
import torch

from .. import CraterDetector


def get_parser():
    parser = argparse.ArgumentParser(
        description="Save model weights to file.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "--run_id",
        type=str,
        default=None,
        nargs="?",
        help="Resume from MLflow run checkpoint",
    )
    parser.add_argument(
        "--output_name",
        type=str,
        default="CraterRCNN",
        nargs="?",
        help="Resume from MLflow run checkpoint",
    )
    return parser


if __name__ == "__main__":
    parser = get_parser()
    args = parser.parse_args()
    mlflow.set_tracking_uri("http://localhost:5000/")
    mlflow.set_experiment("crater-detection")

    model = CraterDetector()
    checkpoint = mlflow.pytorch.load_state_dict(rf"runs:/{args.run_id}/checkpoint")
    model.load_state_dict(checkpoint["model_state_dict"])

    torch.save(model.state_dict(), rf"blobs/{args.output_name}.pth")

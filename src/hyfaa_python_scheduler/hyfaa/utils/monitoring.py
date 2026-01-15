# encoding: utf-8
#
# Utility functions related to the monitoring activity
# (Handling metrics)

from prometheus_client import write_to_textfile, push_to_gateway
from urllib.error import URLError

def write_prometheus_metrics(registry, metrics_filepath=None, pushgateway_url=None):
    """
    Write the metrics generated in the given registry.
    Supports both writing the metrics into a text file (in server's filesystem) or pushing them to a prometheus
    pushgateway. Or both.
    :param registry: the registry containing the metrics
    :param metrics_filepath: path to the textfile. It not provided, no text export will be made.
    :param pushgateway_url: URL of the pushgateway. If not provided, the metrics won't be exported.
    :return:
    """
    if pushgateway_url:
        try:
            push_to_gateway(pushgateway_url, job='preprocessing_forcing', registry=registry)
        except URLError as e:
            print(f"Error while pushing metrics. Cannot connect to pushgateway at {pushgateway_url}. Error message is {e}")

    if metrics_filepath:
        try:
            write_to_textfile(metrics_filepath, registry)
        except IOError as e:
            print(f"Error while writing metrics to file {metrics_filepath}. Error message is {e}")

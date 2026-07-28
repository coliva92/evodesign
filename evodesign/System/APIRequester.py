import time
from typing import Callable, Optional

import requests

from ..RetrievableSettings import RetrievableSettings
from .Exceptions import *




class APIRequester(RetrievableSettings):
    """
    Internal utility for sending HTTP requests to remote APIs.
    """

    def __init__(
        self,
        url: str,
        json_request_key: Optional[str] = None,
        json_response_key: Optional[str] = None,
        sleep_time: float = 1.5,
        connection_timeout: int = 30,
        verify: bool = False,
    ) -> None:
        """
        Args:
            url: The remote API endpoint.
            json_request_key: Optional key to wrap outgoing payloads (required by some backend schemas).
            json_response_key: Optional key to unwrap incoming JSON responses.
            sleep_time: Delay in-between requests to prevent rate-limiting (in seconds).
            connection_timeout: Max time to wait for a connection (in seconds).
            verify: Whether to verify SSL certificates.
        """
        self.url = url
        self.json_request_key = json_request_key
        self.json_response_key = json_response_key
        self.sleep_time = sleep_time
        self.connection_timeout = connection_timeout
        self.verify = verify
        return



    def post(self, payload_data):
        return self._send_request(payload_data, requests.post)



    def get(self, payload_data):
        return self._send_request(payload_data, requests.get)



    def _send_request(
        self,
        payload_data,
        method: Callable,
    ):
        if self.sleep_time > 0.0:
            time.sleep(self.sleep_time)
        if self.json_request_key is not None:
            response = method(
                self.url,
                json={self.json_request_key: payload_data},
                timeout=self.connection_timeout,
                verify=self.verify,
            )
        else:
            response = method(
                self.url,
                data=payload_data,
                timeout=self.connection_timeout,
                verify=self.verify,
            )
        if response.status_code == 400:
            raise HttpBadRequest
        elif response.status_code == 403:
            raise HttpForbidden
        elif response.status_code == 504:
            raise HttpGatewayTimeout
        elif response.status_code == 500:
            raise HttpInternalServerError
        elif response.status_code != 200:
            print(response.status_code)
            print(response.content.decode())
            raise HttpUnknownError
        if self.json_response_key is not None:
            return response.json()[self.json_response_key]
        return response.content.decode()

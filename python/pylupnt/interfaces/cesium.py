from fastapi import FastAPI, Response
from pathlib import Path
import threading
import uvicorn
import pylupnt as pnt
import json
import numpy as np
from datetime import datetime, timedelta
import os


class CesiumViewer:
    def __init__(self, host="127.0.0.1", port=8080):
        from pylupnt.core.base import BASEDIR

        self.static_path = BASEDIR / "local" / "cesium"
        os.makedirs(self.static_path, exist_ok=True)
        self.host = host
        self.port = port
        self.thread = None
        self.server = None
        self.config = None
        self.entities = []  # Store entity data

        self.cesium_token = os.getenv("CESIUM_TOKEN")
        if not self.cesium_token:
            pnt.Logger.warning("CESIUM_TOKEN environment variable is not set", name="Cesium")
            return

        from pylupnt.core.pylupnt_utils import LUPNT_DATA_PATH

        index_template_path = LUPNT_DATA_PATH / "cesium" / "index.html"
        index_content = index_template_path.read_text()
        index_content = index_content.replace("{{ CESIUM_TOKEN }}", self.cesium_token)
        (self.static_path / "index.html").write_text(index_content)

        self.app = FastAPI()
        self._setup_routes()
        self.start()

    def _setup_routes(self):
        @self.app.get("/")
        def serve_index():
            file = self.static_path / "index.html"
            if not file.exists():
                return Response(status_code=404)
            return Response(content=file.read_text(), media_type="text/html")

        @self.app.get("/entities/{entity_id}.czml")
        def get_entity_czml(entity_id: str):
            """Return a single entity as CZML"""
            entity = next((e for e in self.entities if e["id"] == entity_id), None)
            if not entity:
                return Response(status_code=404, content=f"Entity {entity_id} not found")

            czml_data = self._create_czml_data_for_entity(entity)
            return Response(
                content=czml_data,
                media_type="application/json",
                headers={"Content-Type": "application/json"},
            )

        @self.app.get("/entities/all")
        def get_all_entities():
            """Return all entities with their body information"""
            if not self.entities:
                return {"entities": []}

            entities_info = []
            for entity in self.entities:
                # Convert BodyId enum to integer for JSON serialization
                body_id = entity.get("body_id", pnt.BodyId.EARTH)
                body_id_int = int(body_id) if hasattr(body_id, "__int__") else int(body_id.value)
                entities_info.append({"entity_id": entity["id"], "body_id": body_id_int})

            return {"entities": entities_info}

        @self.app.get("/entities/timerange")
        def get_time_range():
            """Return the global time range for all entities"""
            if not self.entities:
                return {"start_time": None, "end_time": None}

            global_start_time = None
            global_end_time = None

            for entity in self.entities:
                start_time = datetime.fromisoformat(
                    entity["initial_time_utc"].replace("Z", "+00:00")
                )
                stop_time = start_time + timedelta(seconds=entity["times"][-1])

                if global_start_time is None or start_time < global_start_time:
                    global_start_time = start_time
                if global_end_time is None or stop_time > global_end_time:
                    global_end_time = stop_time

            return {
                "start_time": (
                    global_start_time.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z"
                    if global_start_time
                    else None
                ),
                "end_time": (
                    global_end_time.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z"
                    if global_end_time
                    else None
                ),
            }

    def add_entity(
        self,
        times,
        positions,
        initial_time_utc,
        frame,
        entity_id="Entity",
        name="Entity",
        color=(255, 0, 255),
        description="Entity orbit",
        body_id=pnt.BodyId.EARTH,
        size=5,
    ):
        """
        Add an entity to the viewer

        Parameters:
        -----------
        times : array-like or float
            Array of times in seconds (for satellites) or single time (for ground entities)
        positions : array-like or list
            Array of positions in meters, shape (n, 3) for n time steps (for satellites)
            or single position [x, y, z] in meters (for ground entities)
        initial_time_utc : str
            ISO string for the initial UTC time (e.g., '2024-06-01T12:00:00Z')
        entity_id : str
            Unique identifier for the entity
        name : str
            Display name for the entity
        color : tuple
            RGB color tuple (r, g, b) for the entity and its path
        description : str
            Description of the entity
        body_id : pnt.BodyId
            Reference body for the entity (e.g., pnt.BodyId.EARTH, pnt.BodyId.MOON)
        size : int
            Size of the entity point
        frame : str
            Reference frame for the entity (e.g., "INERTIAL", "FIXED")
        """
        # Convert to numpy arrays if needed
        times = np.asarray(times)
        positions = np.asarray(positions) * 1e3  # Convert to meters

        if body_id == pnt.BodyId.MOON and frame == "INERTIAL":
            raise ValueError("Moon inertial frame is not supported")

        # Handle single position (ground entity)
        if positions.ndim == 1 and len(positions) == 3:
            # Single position - reshape to (1, 3)
            positions = positions.reshape(1, 3)
            times = np.asarray([times]) if np.isscalar(times) else times.reshape(1)

        # Ensure positions has correct shape
        if positions.ndim == 1:
            positions = positions.reshape(-1, 3)

        if len(times) != len(positions):
            raise ValueError("Times and positions arrays must have the same length")

        # Store entity data
        entity_data = {
            "id": entity_id,
            "name": name,
            "times": times.tolist(),
            "positions": positions.tolist(),
            "color": color,
            "description": description,
            "body_id": body_id,
            "initial_time_utc": initial_time_utc,
            "size": size,
            "frame": frame,
        }

        self.entities.append(entity_data)

        # Log the addition
        pnt.Logger.info(
            f"Added {name} to {body_id.name} viewer",
            name="Cesium",
        )

    def _create_czml_data_for_entity(self, entity):
        """Create CZML data for a single entity"""
        start_time = datetime.fromisoformat(entity["initial_time_utc"].replace("Z", "+00:00"))
        stop_time = start_time + timedelta(seconds=entity["times"][-1])

        # Create position data with proper time format
        position_data = []
        for i, (time, pos) in enumerate(zip(entity["times"], entity["positions"])):
            # Calculate the actual time for this position
            current_time = start_time + timedelta(seconds=float(time))
            # Format as ISO string for CZML (use UTC format without timezone offset)
            time_str = current_time.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z"
            position_data.extend([time_str, pos[0], pos[1], pos[2]])

        # Create CZML entity
        czml_entity = {
            "id": entity["id"],
            "name": entity["name"],
            "description": entity["description"],
            "availability": f"{start_time.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3]}Z/{stop_time.strftime('%Y-%m-%dT%H:%M:%S.%f')[:-3]}Z",
            "label": {
                "fillColor": {"rgba": [255, 255, 255, 255]},
                "font": "13pt Lucida Console",
                "horizontalOrigin": "LEFT",
                "outlineColor": {"rgba": [0, 0, 0, 255]},
                "outlineWidth": 3,
                "pixelOffset": {"cartesian2": [20, 0]},
                "style": "FILL_AND_OUTLINE",
                "text": entity["name"],
            },
            "position": {
                "cartesian": position_data,
                "interpolationAlgorithm": "LAGRANGE",
                "interpolationDegree": 1,
                "referenceFrame": entity["frame"],
            },
            "point": {
                "pixelSize": entity["size"],
                "color": {
                    "rgba": [
                        entity["color"][0],
                        entity["color"][1],
                        entity["color"][2],
                        255,
                    ]
                },
                "outlineColor": {
                    "rgba": [
                        entity["color"][0],
                        entity["color"][1],
                        entity["color"][2],
                        255,
                    ]
                },
                "outlineWidth": 2,
                "show": True,
            },
        }

        # Add path for trajectory visualization (only if we have multiple positions)
        if len(entity["times"]) > 1:
            czml_entity["path"] = {
                "show": True,
                "width": 2,
                "material": {
                    "solidColor": {
                        "color": {
                            "rgba": [
                                entity["color"][0],
                                entity["color"][1],
                                entity["color"][2],
                                128,
                            ]
                        }
                    }
                },
                "resolution": 120,
            }

        # Create a simple CZML document
        czml_document = [
            {"id": "document", "name": "Satellite Visualization", "version": "1.0"},
            czml_entity,
        ]
        czml_json = json.dumps(czml_document, indent=2)
        return czml_json

    def clear_entities(self):
        """Clear all entities from the viewer"""
        self.entities.clear()

    def start(self):
        if self.thread and self.thread.is_alive():
            return

        import socket

        def is_port_in_use(host, port):
            with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
                s.settimeout(0.5)
                return s.connect_ex((host, port)) == 0

        def run():
            port = self.port
            # Find a free port before starting the server
            while is_port_in_use(self.host, port):
                port += 1

            self.config = uvicorn.Config(self.app, host=self.host, port=port, log_level="warning")
            self.server = uvicorn.Server(self.config)
            pnt.Logger.info(
                f"Server running at http://{self.host}:{port}",
                name="Cesium",
            )
            self.server.run()

        self.thread = threading.Thread(target=run, daemon=True)
        self.thread.start()

    def stop(self):
        if self.server and self.server.should_exit is False:
            self.server.should_exit = True
        if self.thread:
            self.thread.join(timeout=2)

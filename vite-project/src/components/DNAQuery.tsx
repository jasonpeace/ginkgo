import { useEffect, useState } from "react";
import { Button } from "react-bootstrap";
import { useInterval } from "useHooks-ts";
import Row from "react-bootstrap/Row";
import Col from "react-bootstrap/Col";
import Table from "react-bootstrap/Table";
import Container from "react-bootstrap/Container";

import "bootstrap/dist/css/bootstrap.min.css";
import "./DNAQuery.css";

type AlignmentRequest = {
  id: number;
  dna_string: string;
  status: string;
  date_submitted: string;
  date_updated: string;
};

function DNAQuery() {
  const [newDNAQuery, setNewDNAQuery] = useState<string>();
  const [alignmentRequests, setAlignmentRequests] = useState([]);

  //Fetch the initial data for the DNA Query Table
  useEffect(() => {
    fetch_all_alignment_requests();
  }, []);

  //This hook will call the fetch function every 30 seconds to update
  //the submitted DNA Query Table
  useInterval(() => {
    fetch_all_alignment_requests();
  }, 30000);

  //Fetches the data for the DNA Query Table. Note there is no pagination or filtering
  function fetch_all_alignment_requests() {
    fetch("/api/requests", {
      method: "GET",
    })
      .then((response) => response.json())
      .then((data) => {
        setAlignmentRequests(data);
      });
  }

  function handleNewDNAQuery(event: React.ChangeEvent<HTMLTextAreaElement>) {
    const data: string = event.target.value;
    const limited_data: string = extract(data, "[actgACTG]+");
    setNewDNAQuery(limited_data);
  }

  const extract = (str: string, pattern: string) =>
    (str.match(pattern) || []).pop() || "";

  //When the user hits submit, this takes the DNA query and sends it to the backend.
  //This will refresh the DNA Query Table.
  const handleSubmit = () => {
    fetch("/api/request", {
      method: "POST",
      headers: {
        Accept: "application/json",
        "Content-Type": "application/json",
      },
      body: JSON.stringify({ dna_string: newDNAQuery }),
    })
      .then(function () {
        setNewDNAQuery("");
        fetch_all_alignment_requests();
      })
      .catch(function (res) {
        console.log(res);
      });
  };

  return (
    <>
      <Container>
        <Row>
          <Col>
            <div className="d-grid gap-2">
              <textarea
                name="DNAQueryForm"
                onChange={handleNewDNAQuery}
                value={newDNAQuery}
                className="dna-query"
              />
            </div>
          </Col>
          <Col>
            <p>Please enter your DNA query in the text field to the left.</p>
            <p>Note: only ACTG in upper or lower case is allowed.</p>
            <Button disabled={!newDNAQuery} onClick={handleSubmit} size="sm">
              Submit DNA
            </Button>
          </Col>
        </Row>
      </Container>

      {alignmentRequests.length > 0 && (
        <div className="dna-table">
          <Table striped bordered hover>
            <thead>
              <tr>
                <td>ID</td>
                <td>DNA String</td>
                <td>Status</td>
                <td>Date Submitted</td>
                <td>Date Updated</td>
                <td>Actions</td>
              </tr>
            </thead>
            <tbody>
              {alignmentRequests.map((request: AlignmentRequest) => {
                return (
                  <tr key={request.id}>
                    <td className="table-id">{request.id}</td>
                    <td className="table-dna">{request.dna_string}</td>
                    <td className="table-status">{request.status}</td>
                    <td className="table-date">{request.date_submitted}</td>
                    <td className="table-date">{request.date_updated}</td>
                    <td className="table-actions">
                      <Button
                        onClick={() => {
                          location.href = "/detail/" + request.id;
                        }}
                      >
                        Details
                      </Button>
                    </td>
                  </tr>
                );
              })}
            </tbody>
          </Table>
        </div>
      )}
    </>
  );
}

export default DNAQuery;

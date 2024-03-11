import { useState, useEffect } from "react";

const RequestDetail = ({ request }) => {
  const [request, setRequest] = useState(request);

  function handleUpdateRequest(data) {
    setRequest(data);
    console.log(data);
  }

  const [result, setResult] = useState();

  function handleUpdateResult(data) {
    setResult(data);
    console.log(data);
  }



  return (
    <>
      <table>
        <thead>
          <tr>
            <td>DNA String</td>
            <td>Status</td>
            <td>Date Submitted</td>
            <td>Date Updated</td>
          </tr>
        </thead>
        <tbody>
          <tr>
            <td>{request.DNAString}</td>
            <td>{request.status}</td>
            <td>{request.dateSubmitted}</td>
            <td>{request.dateUpdated}</td>
          </tr>
        </tbody>
      </table>

      <table>
        <thead>
          <tr>
            <td>Protein id</td>
            <td>Alignment Result</td>
          </tr>
        </thead>
        <tbody>
          <tr>
            <td>{result.proteinId}</td>
            <td>{result.alignmentResult}</td>
          </tr>
        </tbody>
      </table>
    </>
  );
};

export default SlideToggle;
